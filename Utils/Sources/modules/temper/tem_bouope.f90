subroutine tem_bouope()
  !------------------------------------------------------------------------
  !****f* Temper/tem_bouope
  ! NAME 
  !    tem_bouope
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, boundary operations
  ! USES
  ! USED BY
  !    tem_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper 
  use def_domain
  use def_temper
  use mod_ADR,    only : ADR_add_sgs_or_bubble
  use mod_solver, only : solver_assemble_element_matrix
  use mod_matrix, only : matrix_assemble_element_RHS
  implicit none
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: elvel(ndime,mnode),gptem(mgaus)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: bocod(ndime,mnodb),botem(mnodb)
  real(rp)    :: bovel(ndime,mnodb)
  real(rp)    :: bogrc(ndime,mnodb)
  integer(ip) :: ielem,inodb,ipoin,kfl_gobou
  integer(ip) :: igaus,igaub,iboun,pblty,idime
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  integer(ip) :: dummi
  real(rp)    :: eucta,tmatr,gbsur,gpdet,adotn
  real(rp)    :: gbsph(mgaub),gbden(mgaub),gbcon(mgaub),gbvis(mgaub),gbtem,gbvel(3)
  real(rp)    :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus)
  real(rp)    :: gprea(mgaus)
  real(rp)    :: gptur(mgaus)                                 ! Turbulent viscosity
  real(rp)    :: arobi,trobi,qrobi,twall,acvis
  real(rp)    :: para1,para2,para3,para4
  real(rp)    :: xmrhs,xmmat
  real(rp)    :: eledd(mnode),gpcar(ndime,mnode,mgaus)
  real(rp)    :: gpcon(mgaus),dummr(ndime*mnode),gpcod, gpvis(mgaus)
  real(rp)    :: xjaci(9),xjacm(9), kinen
  real(rp)    :: temex,velex(3),veave(3)
  integer(ip)             :: kk, ifiel
  real(rp)                :: xx, temwa,rhocp, beta

  !
  ! Loop over elements  
  !
  !  call tem_set_heat_flux( bvnat_tem(3,1:nboun,1) )
  !
  boundaries: do iboun=1,nboun

     kfl_gobou = 1

     Neumann_or_Robin: if( (                &
          &  kfl_fixbo_tem(iboun) == 2 .or. &
          &  kfl_fixbo_tem(iboun) == 3 .or. &
          &  kfl_fixbo_tem(iboun) == 4 .or. &
          &  kfl_fixbo_tem(iboun) == 5 .or. &
          &  bemol_tem /= 0.0_rp )  .and. kfl_gobou == 1 ) then

        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)

        if( pelty > 0 ) then

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
           botem(1:pnodb)         = therm(lnodb(1:pnodb,iboun),1)
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           if( kfl_fixbo_tem(iboun) == 5 ) then
              bogrc(1:ndime,1:pnodb) = gradc_tem(1:ndime,lnodb(1:pnodb,iboun))
           end if
           !
           ! GPTEM+GPSGS: Temperature at Gauss point
           !
           call gather(&
                1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
                elmar(pelty)%shape,therm,gptem)
           call ADR_add_sgs_or_bubble(&
                ielem,pgaus,elmar(pelty) % shape_bub,ADR_tem,gptem) 
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
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar) 

           call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
           call ker_proper('CONDU','PGAUB',dummi,iboun,gbcon)
           call ker_proper('SPHEA','PGAUB',dummi,iboun,gbsph)
           call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis)

           if (kfl_fixbo_tem(iboun) == 3) &
                call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape,gpcar)
           ! 
           ! Coupling with turbul
           !
           call tem_turbul(&
                ielem,pnode,pgaus,1_ip,pgaus,elmar(pelty)%shape,gpcon,gpsph,gpdif,dummr,gpden,gptur)

           !
           ! Reaction term
           !
           call tem_elmrea( &
                1_ip,pnode,pgaus,1_ip,pgaus,elvel,gpden,gpcar,&
                gprea)

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
              gbsur = elmar(pblty)%weigp(igaub)*eucta 
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
              !
              ! Cylindrical coordinates
              !
              if( kfl_naxis == 1 ) then
                 gpcod = 0.0_rp
                 do inodb = 1,pnodb
                    gpcod = gpcod + bocod(1,inodb) * elmar(pblty) % shape(inodb,igaub)
                 end do
                 gbsur = gbsur * gpcod * twopi
              end if
              !
              ! Robin: a.n
              !
              if( kfl_advec_tem /= 0 .and. bemol_tem /= 0.0_rp ) then
                 if( kfl_advec_tem == 1 ) then
                    do inodb=1,pnodb
                       ipoin=lnodb(inodb,iboun)
                       bovel(1:ndime,inodb)=veloc(1:ndime,ipoin,1)
                    end do
                 else if(kfl_advec_tem>=2) then
                    call tem_velfun(pnodb,bocod,bovel)
                 end if
                 gbvel=0.0_rp
                 do inodb=1,pnodb
                    gbvel(1:ndime)=gbvel(1:ndime)&
                         +bovel(1:ndime,inodb)*elmar(pblty)%shape(inodb,igaub)
                 end do
                 adotn = dot_product(gbvel(1:ndime), baloc(1:ndime, ndime))

              end if

              qrobi = 0.0_rp
              arobi = 0.0_rp
              trobi = 0.0_rp 

              if( kfl_fixbo_tem(iboun) == 2 ) then
                 !
                 ! Robin condition
                 ! k*grad(T).n = qr+ar*(T-Tr)
                 !
!
                 call tem_set_heat_flux( iboun, bvnat_tem(1:3,iboun,1) )
!
                 arobi = bvnat_tem(1,iboun,1)                        ! ar
                 trobi = bvnat_tem(2,iboun,1)                        ! Tr   
                 qrobi = bvnat_tem(3,iboun,1)                        ! qr    

              else if( kfl_fixbo_tem(iboun) == 4 ) then
                 !
                 ! Augmented Robin condition
                 ! k*grad(T).n = P3+P1*T+P2*(exp(T/P4)-1.0)
                 !
                 para1 = bvnat_tem(1,iboun,1)                        ! P1
                 para2 = bvnat_tem(2,iboun,1)                        ! P2   
                 para3 = bvnat_tem(3,iboun,1)                        ! P3    
                 para4 = bvnat_tem(4,iboun,1)                        ! P4 
                 gbtem = 0.0_rp
                 do inodb = 1,pnodb                  
                    gbtem = gbtem + botem(inodb) *elmar(pblty) % shape(inodb,igaub)
                 end do
                 qrobi = -para3-para2*(exp(gbtem/para4)-1.0_rp)
                 arobi = -para1
                 trobi = 0.0_rp 

              else if( kfl_fixbo_tem(iboun) == 3 ) then
                 !
                 ! Law of the wall
                 !              
                 ! k*grad(T).n = -(rho*cp*u/T+)*(T-Tw)
                 !                  

                 ! CRITERIA
                 ! CODES BOUNDARIES
                 !   BOUNDCODE   CODE(3)   VALUE    IFIEL
                 ! END_CODES
                 ! BVNAT(1)= VALUE  
                 ! BVNAT(2) =IFIEL (from which field)


                 ifiel = int(bvnat_tem(2, iboun,1))
                 if (ifiel==0) then ! constant value for field
                    twall = bvnat_tem(1, iboun,1)
                 else  ! ifiel /= 0 
                    twall = 0.0_rp
                    kk = k_tran_fiel(ifiel) !indicates to which interval the current time belongs.
                    xx = x_tran_fiel(ifiel) !indicates the position between the begining and end of the interval. 

                    do inodb =1, pnodb
                       ipoin = lnodb(inodb,iboun)
                       temwa = xfiel(ifiel) % a(1,ipoin,kk) * xx + xfiel(ifiel) % a(1,ipoin,kk+1) * (1.0_rp-xx)
                       twall= twall + temwa * elmar(pblty)%shape(inodb,igaub)
                    end do

                 end if
                 !                    twall =  265.0 - 6.9444e-05*cutim ! gabls1
                 if( kfl_rough > 0 )then
                    rough_dom = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)                
                       rough_dom=rough_dom+ rough(ipoin)*&
                            elmar(pblty)%shape(inodb,igaub)
                    end do
                 end if

                 kinen = 0.0_rp
                 if( kfl_ustar == 2 ) then 
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       kinen = kinen + untur(1,ipoin,1) * elmar(pblty)%shape(inodb,igaub)
                    end do
                    if (kfl_logva==1) kinen = exp(kinen)                    
                 end if

                 if ( kfl_waexl_ker == 1_ip ) then !if exchange location for wall law
                    temex = temel_ker(1,lexlo_ker(igaub,iboun))
                    velex(1:ndime) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))
                 else
                    temex = 0.0_rp
                    velex(1:ndime) = 0.0_rp
                 end if

                 if ( kfl_wlaav_ker == 1_ip ) then ! Time-averaged velocity for wall law
                    veave(1:ndime) = velav_ker(1:ndime,igaub,iboun)
                 else
                    veave(1:ndime) = 0.0_rp
                 end if
                 call tem_bouwal(&
                      lboel(1,iboun),iboun, elmar(pblty)%shape(1,igaub),pnodb,&
                      pnode,gbden(igaub),gbvis(igaub),&
                      gbsph(igaub),gbcon(igaub),prtur_tem,twall,rough_dom,&
                      kinen,baloc,velex,veave,arobi,trobi,qrobi, gbvel)

              else if( kfl_fixbo_tem(iboun) == 5 ) then
                 !
                 ! Robin condition for water vapor model
                 ! 
                 gbtem = 0.0_rp
                 do inodb = 1,pnodb
                    gbtem = gbtem + botem(inodb) * elmar(pblty) % shape(inodb,igaub)
                 end do
                 call tem_watvap(&
                      pnodb,ndime,baloc(1,ndime),bogrc,&
                      elmar(pblty) % shape(1,igaub),&
                      gbtem,qrobi,arobi,trobi)

              elseif( kfl_fixbo_tem(iboun) == 20 ) then
                 !
                 ! Stable outflow condition (implicit)
                 ! beta*rho * {u · n}_{-} gptem
                 !
                 ! boundary velocity

                 beta = bvnat_tem(1,iboun,1)             
                 !
                 ! implicit boundary term to matrix
                 !                 
                 arobi =  min(adotn, 0.0_rp)*gbden(igaub)*gbsph(igaub)*beta

              end if
              rhocp = gbsph(igaub)*gbden(igaub)
              if (  kfl_fixbo_tem(iboun) == 3.and.&
                   kfl_waexl_imp_ker == 0 .and.kfl_waexl_ker==1) then
                 ! EXPLICIT EXCHANGE LOCATION
                 xmrhs =  qrobi - arobi*trobi                         !  qr - ar * Tr
                 xmrhs =  xmrhs + (arobi - bemol_tem*adotn*rhocp)*temex     ! -ar + bemol(u.n)
                 xmmat =  0.0_rp 
              else
                 xmrhs =  qrobi - arobi*trobi                                !  qr - ar * Tr
                 xmmat =- arobi + bemol_tem*adotn*rhocp   ! -ar + bemol(u.n)
              end if
              call tem_boumat(&
                   pnode,pnodb,igaub, iboun,lboel(1,iboun),lelbo(iboun),xmmat,xmrhs,&
                   elmar(pblty)%shape(1,igaub),gbsur,elmat,elrhs)

           end do gauss_points
           !
           ! Prescribe Dirichlet boundary conditions
           !
           if( solve(1) % kfl_iffix == 0 ) &
                call tem_elmdir(&
                pnode,lnods(1,ielem),elmat,elrhs,ielem)
           !
           ! Assembly
           !
           !call assrhs(solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
           call solver_assemble_element_matrix(&
                solve_sol,1_ip,pnode,pnode,ielem,lnods(:,ielem),elmat,amatr)

           !call assmat(&
           !     solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
           !     ielem,lnods(1,ielem),elmat,amatr)

        end if

     end if Neumann_or_Robin

  end do boundaries


contains

  !-----------------------------------------------------------------------||---!
  subroutine tem_set_heat_flux( i_boun, param )

  use def_domain,           only: npoin, nboun
  use mod_ker_space_time_function, only : ker_space_time_function 
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in   ) :: i_boun 
  real(rp),    intent(inout) :: param(:) ! ar*(T-Tr) + qr  
  logical(ip) :: is_robin, is_code
  real(rp)    :: dummy=0.0_rp, newY=huge(0.0_rp)  
  real(rp)    :: hconv, Tref, q0 
  !-----------------------------------------------------------------------||---!
  !! JMZA. 2018MAR25
  !!
  !! kfl_icodb ->  1 
  !! CODES, BOUNDARIES
  !!   1  2  0.0  0.0  0.0
  !! END_CODES_BOUNDARIES
  !!
  !! kfl_conbc_tem -> 0                             ____ kfl_funbo!   
  !! BOUNDARY_CONDITIONS , NONCO                   / 
  !!   IDboundary kfl_fixbo, val1 val2 ..., FUNCT=IDFUNC  
  !!   ...       
  !!
  !print *, "+ kfl_icodb, kfl_conbc_tem, number_space_time_function ",  kfl_icodb>0, kfl_conbc_tem==0, number_space_time_function 
  !print *, " ", shape(kfl_funbo), shape(kfl_fixbo)
    
  if( (kfl_icodb>0).and.(kfl_conbc_tem==0) ) then
    hconv = param(1)
    Tref  = param(2)
    q0    = param(3)
    call ker_space_time_function(-kfl_funbo_tem(i_boun),dummy,dummy,dummy,cutim,newY)

    is_robin = kfl_fixbo_tem(i_boun)==2 
    if( is_robin ) then
      hconv = param(1)
      Tref  = param(2)
      q0    = param(3)
      param = (/ hconv, newY, q0 /) 

      !print *, i_boun, cutim, newY, param(:) 
      !print *, i_boun, kfl_codbo(i_boun), kfl_fixbo_tem(i_boun), kfl_funbo_tem(i_boun) 

    endif 

    
  endif 
!-----------------------------------------------------------------------||---!

end subroutine tem_set_heat_flux
  !-----------------------------------------------------------------------||---!
  ! + tem_reabcs.f90 
  ! |_ domain/mod_opebcs.f90. boundary_conditions_read_boundary_codes  
  !    reacod(200)
  !   |_ domain/reacod.f90. 
  !
  ! + tem_updbcs.f90
  ! |_ kfl_conbc_tem==0 ! Non-constant bc 
  !   |_ kfl_fixbo_tem(iboun) /= 0; kfl_funbo_tem(iboun) > 0   
  !      tenew = bvnat_tem(ipnat,iboun,2) * funcre( funpa_tem(1,kfl_funbo_tem(iboun)),6,kfl_funty_tem(kfl_funbo_tem(iboun)), cutim)
  !   |_ number_space_time_function >0; kfl_funno_tem(ipoin) < 0 
  !      call( -kfl_funno_tem(ipoin),..., tenew ) 
  !      bvess_tem(1,ipoin,1) = tenew    
  !
  ! + tem_inibcs.f90. 
  !   if( kfl_icodb > 0      ) call reacod(20_ip) 
  !   if( kfl_conbc_tem == 0 ) call tem_membcs(2_ip)   
  ! |_ domain/reacod.f90 
  !    if( itask == 20 ) ! Boundary codes  
  !      kfl_funbo(iboun) ->  
  ! |_ tem_membcs.f90 ! Non-constant b.c.'s : Functions 
  !    allocate( kfl_funbo_tem(nboun) )   
  ! |_ 
  !      !! JMZA. 2018MAR25 
  !      if(kfl_conbc_tem==0) then
  !         iffun     =  1
  !         kfl_funbo => kfl_funbo_tem
  !      else
  !         iffun      =  0
  !      end if
  !     
  !-----------------------------------------------------------------------||---!


end subroutine tem_bouope
