subroutine tem_updthe(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_updthe
  ! NAME 
  !    tem_updthe
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    tem_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use def_domain
  use def_solver
  use def_kermod,         only : gasco
  use mod_communications, only : PAR_SUM
  use mod_ADR,            only : TRAPEZOIDAL
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: baloc(ndime,ndime),bovel(ndime,mnodb) , botem(mnodb)
  real(rp)                :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                :: eltem(mnode),gptem,gpvol(mgaus)
  real(rp)                :: gbsur,eucta,gbvel(ndime),gphes(ntens,mnode), gbtem
  real(rp)                :: gpcar(ndime,mnode,mgaus), fact0, fact1,bmflx
  integer(ip)             :: ielem,inode,ipoin,idime, itime, jtime
  integer(ip)             :: pnode,pgaus,iboun,igaub,inodb,igaus
  integer(ip)             :: pelty,pmate,pblty,pnodb,pgaub,pevat, kelem
  real(rp),    save       :: rinit=0.0_rp

  if( abs(itask) == 1 ) then  ! update thrpr
     !
     ! Integrate thermodynamic pressure in time, used for closed systems with inflow/outflow
     !
     if( kfl_prthe == 2 ) then
        !
        ! inverse of temp and surface flux integration 
        !
        bmflx        = 0.0_rp
        vinvt_tem(1) = 0.0_rp
        if( INOTMASTER ) then
           call tem_smflow(bmflx)
           call tem_vinvte(vinvt_tem(1))
        end if
        !
        ! Parall: Sum over all subdomains
        !
        call PAR_SUM(bmflx,         'IN MY ZONE')
        call PAR_SUM(1_ip,vinvt_tem,'IN MY ZONE')

        if( itask == -1 ) then
           ! 
           ! if initialization task
           ! if initial mass was set, then
           !
           if( imass > zeror )   then  
              !
              ! Set by archive
              !
              vinvt_tem(2) = imass * gasco / prthe(4) * vodom              
           else 
              !
              !             initial mass was not fixed,                  
              !
              !              call runend('TEM_UPDTHE:PLEASE, SET THE INITIAL MASS')
              print *, 'TEM_UPDTHE:INITIAL MASS NOT FIXED FOR THPRE!!'
              vinvt_tem(2) = vinvt_tem(1) 
           end if

        end if
        !
        ! Update thermodynamic pressure
        !
        !  +-             +-          
        !  | d rho        |    
        !  | ----- dv +   | rho*u.n ds = 0
        !  |  dt          |         
        ! -+V            -+S          
        !
        !  rho = p/RT
        !
        !       +-             +-          
        !  d    | pth          | u.n   
        !  --   | --- dv + pth | --- ds = 0
        !  dt   |  T           |  T     
        !      -+V           -+S          
        !
        !  d
        !  -- [ pth * vinvt_tem ] + pth * bmflx = 0
        !  dt
        !
        !  [ 1/dt*c_1*vinvt_tem_1 + bmflx ] * pth + sum_i=2 1/dt * ( c_i*pth_i*vinvt_tem_i ) = 0
        !  [ c_1*vinvt_tem_1 + bmflx/dtinv ] * pth + sum_i=2 ( c_i*pth_i*vinvt_tem_i ) = 0
        !  [ c_1*vinvt_tem_1 + bmflx/dtinv ] * pth + fact0 = 0
        !
        fact0 = 0.0_rp
        do itime = 2,ADR_tem % ntime
           fact0 = fact0 -  ADR_tem % time_parameters(itime) * prthe(itime) * vinvt_tem(itime)
        end do
        prthe(1) = fact0 / (vinvt_tem(1)*ADR_tem % time_parameters(1) + bmflx/ADR_tem % dtinv)
        !
        ! Update thermodynamic pressure time derivative
        !
        dpthe = 0.0_rp
        do itime =1,ADR_tem % ntime
           dpthe = dpthe + ADR_tem % time_parameters(itime) * prthe(itime)
        end do
        dpthe     = dpthe * ADR_tem % dtinv
        xmass_tem = prthe(1) * vinvt_tem(1) / gasco      

     else if( kfl_prthe == 1 ) then   ! Thermodynamic pressure from mass conservation, 
        !used in closed systems without flow
        !
        ! inverse of temperature volume integration
        !
        vinvt_tem(1)=0.0_rp           
        if(INOTMASTER) call tem_vinvte(vinvt_tem(1)) 
        call pararr('SUM',0_ip,1_ip,vinvt_tem(1))

        if(itask==-1) then  ! initialization task
           if (imass > zeror)   then   ! if initial mass was set in datafile
              vinvt_tem(2) = imass*gasco/prthe(4)*vodom 
           else
              vinvt_tem(2) = vinvt_tem(1)
           end if        
        end if
        !
        ! Update thermodynamic pressure and its time derivative
        !
        prthe(1)  = prthe(4)*vinvt_tem(2)/vinvt_tem(1)
        
        dpthe = 0.0_rp
        do itime = 1,ADR_tem % ntime
           dpthe = dpthe + ADR_tem % time_parameters(itime) * prthe(itime)
        end do
        dpthe = dpthe * ADR_tem % dtinv
        ! total  mass prhte/R* integr(1/T dvolu)
        xmass_tem = prthe(1)*vinvt_tem(1)/gasco
     end if

  else if( itask == 2 ) then
     !
     ! Crank-Nicolson
     !
     if( kfl_prthe == 0 ) then
        if( ADR_tem % kfl_time_scheme == TRAPEZOIDAL .and. ADR_tem % kfl_time_order == 2 ) then
           prthe(1) = 2.0_rp * prthe(1) - prthe(2)
        end if
     end if

     do itime = ADR_tem % ntime,2,-1
        prthe(itime) = prthe(itime-1)
        if( kfl_prthe == 2 ) vinvt_tem(itime) = vinvt_tem(itime-1)
     end do
     !do jtime = 2,kfl_tiaor_tem + 1 
     !   itime = kfl_tiaor_tem + 3 - jtime
     !   prthe(itime) = prthe(itime-1)
     !   if( kfl_prthe == 2 ) vinvt_tem(itime) = vinvt_tem(itime-1)
     !end do
  end if

end subroutine tem_updthe

subroutine tem_vinvte(vinvt)
  !-----------------------------------------------------------------------
  !****f* Nastin/tem_vinvt
  ! NAME 
  !    tem_vinvt
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    tem_updthe
  !***
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_temper
  use def_domain
  use def_solver
  implicit none
  real(rp), intent(out)   :: vinvt
  real(rp)                :: baloc(ndime,ndime),bovel(ndime,mnodb) 
  real(rp)                :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                :: eltem(mnode),gptem,gpvol(mgaus)
  real(rp)                :: gbsur,eucta,gphes(ntens,mnode)
  real(rp)                :: gpcar(ndime,mnode,mgaus), fact0, fact1
  integer(ip)             :: ielem,inode,ipoin,idime, itime
  integer(ip)             :: pnode,pgaus,iboun,igaub,inodb,igaus
  integer(ip)             :: pelty,pmate,pblty,pnodb,pgaub, kelem
  logical                 :: kfl_tesgs

  !
  ! Initializations
  !
  vinvt=0.0_rp
  !
  ! Loop on elements
  !

  if (associated(ADR_tem % sgs)) then 
     kfl_tesgs = .true.
  else
     kfl_tesgs = .false.
  end if

  elements: do ielem = 1,nelem
     pelty = ltype(ielem)                  ! Element properties and dimensions
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)
     pmate = 1
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        eltem(inode)=tempe(ipoin,1)
        do idime=1,ndime
           elcod(idime,inode)=coord(idime,ipoin)
        end do
     end do
     if(pmate/=-1) then
        call elmcar(&
             pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
             gphes,ielem)

        do igaus = 1,pgaus
           gptem = 0.0_rp
           do inode = 1,pnode
              gptem = gptem + elmar(pelty) % shape(inode,igaus)*eltem(inode)
           end do
           if( kfl_tesgs ) gptem = gptem + ADR_tem % sgs(ielem)%a(1,igaus,1)
           vinvt = vinvt + gpvol(igaus)/gptem
        end do
     end if
  end do elements

end subroutine tem_vinvte

subroutine tem_smflow(bmflx)
  !-----------------------------------------------------------------------
  !****f* Nastin/tem_vinvt
  ! NAME 
  !    tem_vinvt
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    tem_updthe
  !***
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_temper
  use def_domain
  use def_solver
  implicit none
  real(rp), intent(out)   :: bmflx
  real(rp)                :: baloc(ndime,ndime),bovel(ndime,mnodb) , botem(mnodb)
  real(rp)                :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                :: eltem(mnode),gptem,gpvol(mgaus)
  real(rp)                :: gbsur,eucta,gphes(ntens,mnode), gbvel(ndime), gbtem
  real(rp)                :: gpcar(ndime,mnode,mgaus), fact0, fact1
  integer(ip)             :: ielem,inode,ipoin,idime, itime
  integer(ip)             :: pnode,pgaus,iboun,igaub,inodb,igaus
  integer(ip)             :: pelty,pmate,pblty,pnodb,pgaub, kelem  
  !
  ! Initializations
  !
  bmflx=0.0_rp
  !
  ! Loop on elements
  ! 

  !
  ! Mass M= int_S u.n ds
  !
  do iboun=1,nboun
     pblty = ltypb(iboun) 
     pnodb = nnode(pblty)
     ielem = lelbo(iboun)
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaub = ngaus(pblty)
     pgaus = ngaus(pelty)
     pmate = lmate(ielem)

     do inodb=1,pnodb
        ipoin=lnodb(inodb,iboun)
        do idime=1,ndime
           botem(inodb)       = tempe(ipoin,1) 
           bovel(idime,inodb) = veloc(idime,ipoin,1)
           bocod(idime,inodb) = coord(idime,ipoin)
        end do
     end do
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        do idime=1,ndime
           elcod(idime,inode) = coord(idime,ipoin)
        end do
     end do
     do igaub=1,pgaub
        call bouder(&
             pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
             bocod,baloc,eucta)                                   ! and Jacobian
        gbsur=elmar(pblty)%weigp(igaub)*eucta 
        call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
        do idime=1,ndime
           gbvel(idime)=0.0_rp
        end do
        gbtem =0.0_rp
        do inodb=1,pnodb
           fact0 = elmar(pblty)%shape(inodb,igaub) 
           do idime=1,ndime
              gbvel(idime)=gbvel(idime)&
                   +bovel(idime,inodb)*fact0                           
           end do
           gbtem = gbtem + botem(inodb)*fact0
        end do
        do idime=1,ndime
           bmflx = bmflx +gbvel(idime)*baloc(idime,ndime)*gbsur/gbtem
        end do
     end do
  end do


end subroutine tem_smflow


