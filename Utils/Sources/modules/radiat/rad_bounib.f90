subroutine rad_bounib()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_bounib
  ! NAME 
  !    rad_bounib
  ! DESCRIPTION
  !    ORDER=1:
  !      Radiation Heat Transfer equation, boundary operations
  ! USES
  ! USED BY
  !    rad_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip) :: ielem,inode,ipoin,igaub,inodb,iboun
  integer(ip) :: jnode,idime,dummi
  integer(ip) :: pnodb,pblty,pgaub,pelty,pnode,pmate,porde
  real(rp)    :: baloc(ndime,ndime),tragl(9)
  real(rp)    :: bocod(ndime,mnoib)
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: elvel(ndime,mnode),elcod(ndime,mnode)
  real(rp)    :: eledd(mnode),gpcar(ndime,mnode)
  real(rp)    :: gbcod(3),eucta,gbsur,hleng(3)
  real(rp)    :: deriv(3,64),coloc(3),shaib(mnode)
  real(rp)    :: gpeps,h,adv,rea,dif,gpvel(3),gptem
  real(rp)    :: gpvno,gpsph,gpden,gptau,gpvol,gphes
  real(rp)    :: gprea,gpdif,gpcon,chale(2),chave(3),gpdet
  real(rp)    :: xjaci(9),xjacm(9),gpcib(ndime,mnode)

!!$  if( kfl_immbo==1 .and. INOTMASTER ) then
!!$     !
!!$     ! Loop over elements  
!!$     !
!!$     do iboun=1,nboib
!!$
!!$        pblty=ltyib(iboun)           ! Dimensions
!!$        pnodb=nnode(pblty)
!!$        pgaub=ngaib(pblty)
!!$
!!$        do inodb=1,pnodb             ! Gather 
!!$           ipoin=lnoib(inodb,iboun)
!!$           do idime=1,ndime
!!$              bocod(idime,inodb) = cooib(idime,ipoin)
!!$           end do
!!$        end do
!!$
!!$        do igaub=1,pgaub             
!!$           do idime=1,ndime                        ! GBCOD: Coordinates of IB Gauss points
!!$              gbcod(idime)=0.0_rp
!!$           end do
!!$           do inodb=1,pnodb
!!$              do idime=1,ndime
!!$                 gbcod(idime)=gbcod(idime)+bocod(idime,inodb)&
!!$                      &       *elmar(pblty)%shaib(inodb,igaub)
!!$              end do
!!$           end do
!!$           call bouder(&                           ! GBSUR: Element surface
!!$                pnodb,ndime,ndimb,elmar(pblty)%derib(1,1,igaub),&
!!$                bocod,baloc,eucta)
!!$           gbsur=elmar(pblty)%weiib(igaub)*eucta 
!!$           call chenor(pnode,baloc,bocod,elcod)    ! BALOC: Check normal
!!$           !
!!$           ! IELEM: host element
!!$           !
!!$           call elsest(&
!!$                2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),pelpo,&
!!$                lelpo,lnods,ltype,ltopo,coord,gbcod,relse,&
!!$                ielem,shaib,deriv,coloc,dummi) 
!!$           if(ielem==-1) then
!!$              call runend('IB METHOD: GAUSS POINT OUT OF SUBDOMAIN BOUNDAING BOX')
!!$           else if(ielem==0) then
!!$              call runend('IB METHOD: COULD NOT FIND HOST ELEMENT')
!!$           end if
!!$           ! 
!!$           ! GPEPS: Find epsilon
!!$           !
!!$           pelty = ltype(ielem)
!!$           pnode = nnode(pelty)
!!$           porde = lorde(pelty)
!!$           pmate = 1
!!$           if(nmate>1) pmate=lmate(ielem)
!!$           do inode=1,pnode                                         ! ELCOD 
!!$              ipoin=lnods(inode,ielem)
!!$              do idime=1,ndime
!!$                 elcod(idime,inode) = coord(idime,ipoin)
!!$              end do
!!$           end do
!!$           call elmder(&
!!$                pnode,ndime,deriv,&                                 ! Cartesian derivative
!!$                elcod,gpcib,gpdet,xjacm,xjaci)                      ! and Jacobian
!!$
!!$           call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,& ! HLENG and TRAGL at center of gravity
!!$                hnatu(pelty),hleng)
!!$           call elmchl(tragl,hleng,elcod,elvel,chave,chale,pnode,&  ! CHALE: characteristic length 
!!$                porde,hnatu(pelty),0_ip,kfl_ellen_rad)
!!$           gpvno = 0.0_rp
!!$           call elmcar(&
!!$                pnode,1_ip,0_ip,elmar(pelty)%weicg,elmar(pelty)%shacg,&
!!$                elmar(pelty)%dercg,elmar(pelty)%hescg,elcod,gpvol,gpcar,&
!!$                gphes,ielem)
!!$           call rad_elmpro(& 
!!$                ielem,pmate,pnode,1_ip,1_ip,1_ip,&
!!$                elmar(pelty)%shacg,gpcar,gpdif,gprea,dummr,dummr)
!!$           adv = gpden*gpsph*gpvno
!!$           dif = gpdif
!!$           rea = gprea
!!$           h   =  max(chale(1),chale(2))
!!$           call tauadr(&
!!$                kfl_taust_rad,staco_rad,adv,dif,rea,&
!!$                chale(1),chale(2),gptau)
!!$           gpeps = 1.0_rp/(h*(dtinv_rad+1.0_rp/gptau))
!!$           gpeps = 1e6*gpeps
!!$           !
!!$           ! ELMAT and ELRHS: Matrix and RHS
!!$           !
!!$           pelty=ltype(ielem)
!!$           pnode=nnode(pelty)
!!$           do inode=1,pnode
!!$              elrhs(inode)=0.0_rp
!!$              do jnode=1,pnode
!!$                 elmat(inode,jnode)=0.0_rp
!!$              end do
!!$           end do
!!$           !do inode=1,pnode
!!$           !   elrhs(inode)=elrhs(inode)+gpeps*gbsur*shaib(inode)*0.0_rp
!!$           !   do jnode=1,pnode
!!$           !      elmat(inode,jnode)=elmat(inode,jnode)&
!!$           !           +gpeps*gbsur*shaib(inode)*shaib(jnode)
!!$           !   end do
!!$           !end do
!!$           do inode=1,pnode
!!$              do jnode=1,pnode
!!$                 elmat(inode,jnode)=elmat(inode,jnode)&
!!$                      +gpeps*gbsur*shaib(inode)*shaib(jnode)
!!$                 do idime=1,ndime
!!$                    !elmat(inode,jnode)=elmat(inode,jnode)&
!!$                    !     +1.0_rp*gbsur*shaib(inode)*gpcib(idime,jnode)*baloc(idime,ndime)
!!$                 end do
!!$              end do
!!$           end do
!!$           !
!!$           ! Prescribe Dirichlet boundary conditions
!!$           !
!!$           !call rad_elmdir(&
!!$           !     pnode,lnods(1,ielem),elmat,elrhs)
!!$           !
!!$           ! Assembly
!!$           !
!!$           call assrhs(solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
!!$           call assmat(&
!!$                solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
!!$                lnods(1,ielem),elmat,amatr)
!!$        end do
!!$
!!$     end do
!!$
!!$  end if

end subroutine rad_bounib
