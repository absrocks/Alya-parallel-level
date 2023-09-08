subroutine wav_elmope(order)
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_elmope
  ! NAME 
  !    wav_elmope
  ! DESCRIPTION
  !
  !    Assemble the wave-equation:
  !    
  !    d^2u           +-  1         -+
  !    ---- - kap*div |  --- grad(u) | = kap*f(x)
  !    dt^2           +- rho        -+
  !
  !    ORDER=1 ... Temperature equation, elemental operations:
  !                1. Compute elemental matrix and RHS 
  !                2. Impose Dirichlet boundary conditions
  !                3. Assemble them
  !    ORDER=2 ... Update the subgrid scale
  ! USES
  ! USED BY
  !    wav_matrix
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_wavequ
  implicit none

  integer(ip), intent(in) :: order                          ! 
  integer(ip)             :: ielem,igaus,idime              ! Indices and dimensions
  integer(ip)             :: inode,jnode,ipoin
  integer(ip)             :: pelty,pmate,pnode
  integer(ip)             :: pgaus,plapl,porde
  real(rp)                :: wmat1(ndime,ndime,mnode*mlapl) ! Created here to avoid using an 
  real(rp)                :: d2sdx(ndime,ndime,ndime*mlapl) ! that is called inside do igaus 
  real(rp)                :: xjaci(ndime,ndime) 
  real(rp)                :: xjacm(ndime,ndime) 
  real(rp)                :: elmat(mnode,mnode)             ! Element matrices
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: gpcar(ndime,mnode)             ! dNk/dxj
  real(rp)                :: gphes(ntens,mnode)             ! dNk/dxidxj
  real(rp)                :: gpcod(ndime)                   ! x
  real(rp)                :: gpvol,gpdet,gplap              ! |J|*w,|J|
  real(rp)                :: gpden,gpkap,gpwav(ncomp_wav)   ! rho, kappa, u
  real(rp)                :: gpsou                          ! Q
  real(rp)                :: dtn0,dtn1,gamm0,gamm1,gamm2    ! Time steps
  real(rp)                :: tragl(ndime,ndime)             ! Stabilization
  real(rp)                :: chale(2),dummr(5),tausg,taupr
  real(rp)                :: hleng(3),resid,dtn2
  real(rp)                :: fact1,fact2,gpvel,gpacc

  select case(order)

  case(1)
     !
     ! Left-hand side
     !
     dtn0  = 1.0_rp/dtinv_wav
     dtn1  = dtold(1)
     dtn2  = dtn0*dtn0
     gamm0 = 2.0_rp/(dtn0*(dtn0+dtn1))
     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pgaus=ngaus(pelty)
        porde=lorde(pelty)
        pmate=lmate_wav(ielem)
        do inode=1,pnode
           do jnode=1,pnode
              elmat(inode,jnode)=0.0_rp
           end do
        end do
        !
        ! Gather operations
        !
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
        do igaus=1,pgaus      
           !
           ! Cartesian derivatives, and volume: GPCAR, PGVOL
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           gpvol=elmar(pelty)%weigp(igaus)*gpdet
           ! 
           ! Properties GPDEN and GPKAP
           !
           gpden=densi_wav(1,pmate)
           gpkap=kappa_wav(1,pmate)
           !
           ! Diffusion matrix: 
           ! Newmark:   beta*c^2*grad(Ni).grad(Nj)
           ! Leap-frog: c^2*grad(Ni).grad(Nj)
           !
           fact1=gpkap/gpden*gpvol 
           if(kfl_tisch_wav==3) fact1=fact1*nebet_wav 
           do inode=1,pnode
              do jnode=1,pnode
                 do idime=1,ndime
                    elmat(inode,jnode)=elmat(inode,jnode)&
                         +fact1*gpcar(idime,inode)*gpcar(idime,jnode)
                 end do
              end do
           end do
           if(kfl_tisch_wav==3.and.kfl_massm_wav==1) then
              !
              ! Time term of Newmark: 1/dt^2*Ni*Nj
              !
              fact1=gpvol/dtn2
              do inode=1,pnode
                 fact2=fact1*elmar(pelty)%shape(inode,igaus)
                 do jnode=1,pnode
                    elmat(inode,jnode)=elmat(inode,jnode)&
                         +fact2*elmar(pelty)%shape(jnode,igaus)
                 end do
              end do              
           end if
           if(kfl_timet_wav==2.and.kfl_tisch_wav==2) then 
              !
              ! Time term of implicit Leap-frog
              !
              fact1=gamm0*gpvol
              do inode=1,pnode
                 fact2=fact1*elmar(pelty)%shape(inode,igaus)
                 do jnode=1,pnode
                    elmat(inode,jnode)=elmat(inode,jnode)&
                         +fact1*elmar(pelty)%shape(jnode,igaus)
                 end do
              end do
           end if
        end do
        if(kfl_onnod_wav==1.and.kfl_timet_wav==2) then
           !
           ! Implicit: boundary conditions
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              if(kfl_fixno_wav(ipoin)==1) then
                 do jnode=1,pnode
                    elmat(inode,jnode)=0.0_rp
                 end do
                 elmat(inode,inode)=vmass(ipoin)*dtinv_wav*dtinv_wav
              end if
           end do
        end if
        call assmat(&
             solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
             ielem,lnods(1,ielem),elmat,amatr)
        
     end do

  case(2)
     !
     ! Right-hand side
     !
     dtn0  =  1.0_rp/dtinv_wav
     dtn1  =  dtold(1)
     dtn2  =  dtn0*dtn0
     gamm0 =  2.0_rp/(dtn0*(dtn0+dtn1))
     gamm1 = -2.0_rp/(dtn0*dtn1)
     gamm2 =  2.0_rp/(dtn1*(dtn0+dtn1))
     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pgaus=ngaus(pelty)
        pmate=lmate_wav(ielem)
        !
        ! Gather operations: ELCOD
        !
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
        ! 
        ! Properties GPDEN and GPKAP
        !
        gpden=densi_wav(1,pmate)
        gpkap=kappa_wav(1,pmate)
        !
        ! HLENG and TRAGL at center of gravity
        !
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE
        !
        call elmchl(tragl,hleng,dummr,dummr,dummr,chale,pnode,&
             porde,hnatu(pelty),0_ip,0_ip)

        do igaus=1,pgaus      
           !
           ! Cartesian derivatives, and volume: GPCAR, PGVOL
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           gpvol=elmar(pelty)%weigp(igaus)*gpdet 
           !
           ! Gauss point values
           !       
           gpcod = 0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 gpcod(idime)=gpcod(idime)&
                      +elmar(pelty)%shape(inode,igaus)*elcod(idime,inode) 
              end do
           end do
           !
           ! Source term: kap*f(x)
           !
           call wav_elmsou(ndime,kfl_sourc_wav,sourc_wav,gpcod,oltim,gpsou)
           fact1=gpkap*gpsou*gpvol
           if(kfl_tisch_wav==3) fact1=fact1*nebet_wav
           do inode=1,pnode
              ipoin=lnods(inode,ielem)     
              rhsid(ipoin)=rhsid(ipoin)+elmar(pelty)%shape(inode,igaus)*fact1
           end do
           !
           ! Newmark with consistent matrix
           !
           if(kfl_tisch_wav==3.and.kfl_massm_wav==1) then
              gpwav(1) = 0.0_rp
              gpvel    = 0.0_rp
              gpacc    = 0.0_rp
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 gpwav(1) = gpwav(1) + elmar(pelty)%shape(inode,igaus)*wavam(ipoin,3)
                 gpvel    = gpvel    + elmar(pelty)%shape(inode,igaus)*wavve_wav(ipoin,2)
                 gpacc    = gpacc    + elmar(pelty)%shape(inode,igaus)*wavac_wav(ipoin,2)
              end do
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)     
                 rhsid(ipoin)=rhsid(ipoin)+gpvol*elmar(pelty)%shape(inode,igaus)&
                      *(gpwav(1)/dtn2+dtinv_wav*gpvel+(0.5_rp-nebet_wav)*gpacc)
              end do
           end if
           !
           ! Subgrid scale effects
           !
           if(kfl_subgs_wav==1) then
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)     
                 rhsid(ipoin)=rhsid(ipoin)&
                      -gpvol*elmar(pelty)%shape(inode,igaus)*&
                      (gamm0*wasgs(ielem)%a(igaus,1)&
                      +gamm1*wasgs(ielem)%a(igaus,2)&
                      +gamm2*wasgs(ielem)%a(igaus,3))
              end do
           end if
           !
           ! Assemble time terms
           !
           if(kfl_onnod_wav==1.and.kfl_timet_wav==2.and.kfl_tisch_wav==2) then
              gpwav(1) = 0.0_rp
              gpwav(2) = 0.0_rp
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 gpwav(1)=gpwav(1)+elmar(pelty)%shape(inode,igaus)*wavam(ipoin,3)
                 gpwav(2)=gpwav(2)+elmar(pelty)%shape(inode,igaus)*wavam(ipoin,4)
              end do
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)              
                 if(kfl_fixno_wav(ipoin)==1) then
                    rhsid(ipoin)=rhsid(ipoin)+bvess_wav(ipoin)&
                         *vmass(ipoin)*dtinv_wav*dtinv_wav
                 else
                    rhsid(ipoin)=rhsid(ipoin)-(gpwav(1)*gamm1+gpwav(2)*gamm2)&
                         *elmar(pelty)%shape(inode,igaus)*gpvol
                 end if
              end do
           end if
        end do
     end do

  case(3)
     !
     ! Update subgrid scale WASGS:
     !
     ! g0*u'^{n+1} + g1*u'^{n} + g2*u'^{n-1} + tau^{-1}*u' = R(uh^n)
     ! R(uh^{n+1}) = kap*f(x) - [g0*uh^{n+1} + g1*uh^{n} + g2*uh^{n-1}]-kap*div[1/rho*grad(uh^{n})]
     ! =>
     ! u'^{n} = 1/[g1+tau^{-1}]*[ R(uh^n) - (g0*u'^{n+1} + g2*u'^{n-1}) ]
     !
     dtn0  =  1.0_rp/dtinv_wav
     dtn1  =  dtold(1)
     gamm0 =  2.0_rp/(dtn0*(dtn0+dtn1))
     gamm1 = -2.0_rp/(dtn0*dtn1)
     gamm2 =  2.0_rp/(dtn1*(dtn0+dtn1))
     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pgaus=ngaus(pelty)
        pmate=lmate_wav(ielem)
        plapl=llapl(pelty)
        !
        ! Gather operations
        !
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
        !
        ! HLENG and TRAGL at center of gravity
        !
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE
        !
        call elmchl(tragl,hleng,dummr,dummr,dummr,chale,pnode,&
             porde,hnatu(pelty),0_ip,0_ip)
        !
        ! Properties GPDEN and GPKAP
        !
        gpden = densi_wav(1,pmate)
        gpkap = kappa_wav(1,pmate)

        do igaus=1,pgaus      
           !
           ! Cartesian derivatives, and volume: GPCAR, PGVOL
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           if(plapl==1) &
                call elmhes(&
                elmar(pelty)%heslo(1,1,igaus),gphes,ndime,pnode,ntens,&
                xjaci,d2sdx,elmar(pelty)%deriv(1,1,igaus),&
                elcod)
           gpvol=elmar(pelty)%weigp(igaus)*gpdet
           !
           ! Gauss point values
           !
           gpwav(1) = 0.0_rp
           gpwav(2) = 0.0_rp
           gpwav(3) = 0.0_rp
           do idime=1,ndime
              gpcod(idime)=0.0_rp
           end do           
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              gpwav(1)=gpwav(1)+elmar(pelty)%shape(inode,igaus)*unkno(ipoin)
              gpwav(2)=gpwav(2)+elmar(pelty)%shape(inode,igaus)*wavam(ipoin,3)
              gpwav(3)=gpwav(3)+elmar(pelty)%shape(inode,igaus)*wavam(ipoin,4)
              do idime=1,ndime
                 gpcod(idime)=gpcod(idime)&
                      +elmar(pelty)%shape(inode,igaus)*elcod(idime,inode) 
              end do
           end do
           !
           ! kap*div[1/rho*grad(u)]
           !
           gplap = 0.0_rp
           if(plapl==1) then
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    gplap=gplap+gphes(idime,inode)*unkno(ipoin)
                 end do
              end do
              gplap=gpkap/gpden*gplap
           end if
           !
           ! Source term GPSOU
           !
           call wav_elmsou(ndime,kfl_sourc_wav,sourc_wav,gpcod,oltim,gpsou)
           !
           ! Subgrid scale
           !
           if(1==1) then
              tausg = chale(1)*chale(1)/(12.0_rp*gpkap/gpden)
              resid = gpkap*gpsou-(gamm0*gpwav(1)+gamm1*gpwav(2)+gamm2*gpwav(3))+gplap
              taupr = 1.0_rp/(gamm0+1.0_rp/tausg)
              wasgs(ielem)%a(igaus,1)=&
                   taupr*(resid&
                   -gamm1*wasgs(ielem)%a(igaus,2)&
                   -gamm2*wasgs(ielem)%a(igaus,3))
           else
              tausg = chale(1)*chale(1)/(12.0_rp*gpkap/gpden)
              resid = gpkap*gpsou-(gamm0*gpwav(1)+gamm1*gpwav(2)+gamm2*gpwav(3))+gplap
              wasgs(ielem)%a(igaus,1)=&
                   1.0_rp/gamm0*(resid&
                   -1.0_rp/tausg*wasgs(ielem)%a(igaus,2)&
                   -gamm1*wasgs(ielem)%a(igaus,2)&
                   -gamm2*wasgs(ielem)%a(igaus,3))              
           end if
        end do
     end do

  end select

end subroutine wav_elmope
