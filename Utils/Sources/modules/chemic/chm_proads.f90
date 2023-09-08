subroutine chm_proads()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_proads
  ! NAME 
  !    chm_proads
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    chm_inivar
  !*** 
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_chemic
  use def_solver
  use mod_gradie
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: ielem,igaus,idime,inode,ipoin,pnode,pgaus,pelty,imeth
  integer(ip) :: iiter,niter,izmat,jzmat,jpoin
  real(rp)    :: elmal(mnode,mnode),elrhs(mnode),elcod(ndime,mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus),gpvol(mgaus),gpdet,fact1
  real(rp)    :: xjaci(9),xjacm(9),fact2,eps,xnorm


  call livinf(79_ip,'SOLVE PROTEIN ADSORPTION',modul)

  if( INOTMASTER ) then
     call inisol()
     do ipoin = 1,npoin
        unkno(ipoin) = 0.0_rp
     end do
     do ielem = 1,nelem
        
        if( lelch(ielem) /= ELHOL ) then
           !
           ! Element properties and dimensions
           !
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather operations: ELCOD
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           !
           ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
           !
           do igaus = 1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpdet  
           end do
           !
           ! Compute element matrix ELMAL and assemble LAPLA_TUR
           !
           call chm_elmlap(&
                one,pnode,pgaus,lnods(1,ielem),lelch(ielem),gpcar,&
                elmar(pelty)%shape,gpvol,elmal,elrhs)
           call assmat(&
                1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
                ielem,lnods(1,ielem),elmal,amatr)
           call assrhs(&
                1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
        end if
     end do
  end if
  !
  ! Solve system: Lapl(g)=0, with g = ywalp on nodes with kfl_fixno_tur == 3 or 4
  !
  call solver(rhsid,unkno,amatr,pmatr) 
  !
  ! Compute wall distance: d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !
  if( INOTMASTER ) then
     call memgen(zero,ndime,npoin)
     call gradie(unkno,gevec)
     do ipoin = 1,npoin
        fact1 = 0.0_rp 
        do idime = 1,ndime
           fact1 = fact1 + gevec(idime,ipoin) * gevec(idime,ipoin)
        end do
        fact2 = fact1 + 2.0_rp*max(unkno(ipoin),0.0_rp)
        if( fact2 < 0.0_rp ) then
           call runend('WRONG DISTANCE TO THE WALL')
        else
           unkno(ipoin) = sqrt(fact2) - sqrt(fact1) + unkno(ipoin)
        end if
     end do
     call memgen(two,ndime,npoin)
  end if
  !
  ! Compute p
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        !proad_chm(ipoin) = unkno(ipoin)
        if( unkno(ipoin) > 0.2_rp ) then
           proad_chm(ipoin) = 0.0_rp
        else 
           proad_chm(ipoin) = 0.5_rp/0.2_rp * ( 0.2_rp - unkno(ipoin) )
           proad_chm(ipoin) = min(0.5_rp,proad_chm(ipoin))
           proad_chm(ipoin) = max(0.0_rp,proad_chm(ipoin))
        end if
     end do
  end if

end subroutine chm_proads

subroutine chm_elmlap(&
     itask,pnode,pgaus,lnods,lelch,gpcar,gpsha,gpvol,elmal,elrhs)
  !----------------------------------------------------------------------
  !****f* Chemic/chm_elmlap
  ! NAME 
  !    chm_elmlap
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  !    Itask == 1 the original version for system: Lapl(f)=-1, with f=0 on wall
  !    Itask == 2 use to extend the wall distance from the wall using system: Lapl(g)=0, with g = ywalp on wall
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_chemic, only     :  kfl_fixno_chm
  use def_elmtyp, only     :  ELEXT
  use def_domain, only     :  mnode,nbopo,ndime,coord,lpoty
  implicit none
  integer(ip), intent(in)  :: itask,pnode,pgaus
  integer(ip), intent(in)  :: lelch
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpvol(pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode,kdime,igaus,ipoin
  real(rp)                 :: fact1,bvess,dummr
  !
  ! Initialization
  !
  do jnode=1,pnode
     elrhs(jnode)=0.0_rp
     do inode=1,pnode
        elmal(inode,jnode)=0.0_rp
     end do
  end do
  !
  ! Laplacian matrix: ( grad p , grad q ), and rhs
  !
  do igaus=1,pgaus
     do inode=1,pnode
        do jnode=inode+1,pnode
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gpcar(kdime,inode,igaus)&
                   &     *gpcar(kdime,jnode,igaus)
           end do
           fact1=fact1*gpvol(igaus)
           elmal(inode,jnode)=elmal(inode,jnode)+fact1
           elmal(jnode,inode)=elmal(jnode,inode)+fact1
        end do
        fact1=0.0_rp
        do kdime=1,ndime
           fact1=fact1+gpcar(kdime,inode,igaus)&
                &     *gpcar(kdime,inode,igaus)
        end do
        fact1=fact1*gpvol(igaus)
        elmal(inode,inode)=elmal(inode,inode)+fact1
        elrhs(inode)=elrhs(inode)+gpvol(igaus)*gpsha(inode,igaus)
     end do
  end do
  !
  ! Extension elements
  !
  if( lelch == ELEXT ) then
     call elmext(&
          4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmal,&
          dummr,elrhs,dummr)
  end if
  !
  ! Prescribe Laplacian on walls
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     !if(lpoty(ipoin)/=0.and.kfl_fixno_chm(1,ipoin)==0.and.(coord(2,ipoin)>=0.0_rp.and.coord(2,ipoin)<=1.0_rp)) then
     if(lpoty(ipoin)/=0.and.kfl_fixno_chm(1,ipoin)==0) then
        fact1=elmal(inode,inode)
        do jnode=1,pnode
           elmal(jnode,inode)=0.0_rp              
           elmal(inode,jnode)=0.0_rp              
        end do
        elmal(inode,inode)=fact1
        elrhs(inode)=0.0_rp
     end if
  end do

end subroutine chm_elmlap
